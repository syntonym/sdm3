import venv
import os
import os.path
import sys

join = os.path.join

if not "--virtualenv" in sys.argv:
    # launch virtualenv python
    os.environ["PATH"] = os.path.realpath(join(os.path.curdir, "env", "bin")) + os.pathsep + os.environ["PATH"]
    os.execvp(os.path.basename(sys.executable), [os.path.basename(sys.executable), os.path.realpath(__file__), "--virtualenv"])
else:
    from jupyter_core.command import main
    import os.path
    import re
    path_to_notebook = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Figure.ipynb")
    sys.argv = [re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0]), "notebook", path_to_notebook]
    sys.exit(main())
