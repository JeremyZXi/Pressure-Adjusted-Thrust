"""Read structural mass properties using the installed OpenRocket calculation engine."""
from pathlib import Path
import math
import os
import shutil
import subprocess
import tempfile
import warnings

_JAVA_SOURCE = 'import java.io.File;\nimport java.util.Locale;\nimport info.openrocket.core.startup.OpenRocketCore;\nimport info.openrocket.core.file.GeneralRocketLoader;\nimport info.openrocket.core.masscalc.MassCalculator;\npublic class OrkMass {\n  public static void main(String[] args) throws Exception {\n    OpenRocketCore.initialize();\n    var loader = new GeneralRocketLoader(new File(args[0]));\n    var document = loader.load();\n    var config = document.getRocket().getSelectedConfiguration();\n    var body = MassCalculator.calculateStructure(config);\n    for (var warning : loader.getWarnings()) System.err.println("ORK_WARNING " + warning);\n    System.out.printf(Locale.ROOT, "ORK_MASS %.17g %.17g %.17g %.17g%n",\n      body.getMass(), body.getCM().x, body.getLongitudinalInertia(), body.getRotationalInertia());\n    System.exit(0);\n  }\n}\n'


def load_ork_mass_data(ork_path, *, openrocket_jar=None, java_path=None):
    """Return mass_kg, cg_m, longitudinal_inertia_kg_m2, rotational_inertia_kg_m2.

    Uses the file's selected flight configuration and OpenRocket mass/CG overrides.
    CG is measured from the nose. Inertias are about the structural CG.
    Excludes assigned motors and their propellant. Include liquid-engine hardware
    as structural components, and exclude liquid propellant from those components
    and mass overrides before using this mass with the liquid engine simulation.

    Requires OpenRocket 24.12 and a Java JDK (17+), not just a JRE.
    On macOS the installed OpenRocket and JetBrains runtimes are auto-detected.
    Else pass paths or set OPENROCKET_JAR and OPENROCKET_JAVA.
    """
    path = Path(ork_path).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f'Not an OpenRocket file: {path}')
    jar_value = openrocket_jar or os.environ.get('OPENROCKET_JAR')
    if jar_value:
        jar = Path(jar_value).expanduser().resolve(strict=True)
    else:
        jars = sorted(Path('/Applications/OpenRocket.app/Contents/Resources/app/jar').glob('OpenRocket-*.jar'))
        if len(jars) != 1:
            raise RuntimeError('Set OPENROCKET_JAR to your OpenRocket 24.12 JAR file.')
        jar = jars[0]
    java_value = java_path or os.environ.get('OPENROCKET_JAVA')
    candidates = [str(java_value)] if java_value else [
        str(Path(os.environ['JAVA_HOME']) / 'bin/java') if os.environ.get('JAVA_HOME') else '',
        '/Applications/IntelliJ IDEA.app/Contents/jbr/Contents/Home/bin/java',
        '/Applications/PyCharm.app/Contents/jbr/Contents/Home/bin/java',
        shutil.which('java') or '',
    ]
    java = None
    for candidate in candidates:
        if not candidate:
            continue
        try:
            modules = subprocess.run([candidate, '--list-modules'], capture_output=True, text=True, timeout=15)
        except (OSError, subprocess.TimeoutExpired):
            continue
        if modules.returncode == 0 and 'jdk.compiler@' in modules.stdout:
            java = candidate
            break
    if java is None:
        raise RuntimeError('Set OPENROCKET_JAVA to a JDK 17+ bin/java executable (source compiler required).')
    with tempfile.TemporaryDirectory(prefix='ork-mass-') as directory:
        source = Path(directory) / 'OrkMass.java'
        source.write_text(_JAVA_SOURCE)
        result = subprocess.run(
            [java, '-Djava.awt.headless=true', '-cp', str(jar), str(source), str(path)],
            capture_output=True, text=True, timeout=120,
        )
    if result.returncode:
        raise RuntimeError('OpenRocket mass calculation failed:\n' + result.stderr + result.stdout[-4000:])
    for line in result.stdout.splitlines():
        if line.startswith('ORK_MASS '):
            values = [float(value) for value in line.split()[1:]]
            if len(values) != 4 or not all(math.isfinite(v) for v in values) or values[0] <= 0:
                raise ValueError('OpenRocket returned invalid or zero structural mass properties.')
            for warning in result.stderr.splitlines():
                if warning.startswith('ORK_WARNING '):
                    warnings.warn(warning.removeprefix('ORK_WARNING '), stacklevel=2)
            return dict(zip(('mass_kg', 'cg_m', 'longitudinal_inertia_kg_m2', 'rotational_inertia_kg_m2'), values))
    raise RuntimeError('OpenRocket did not return mass properties.')
