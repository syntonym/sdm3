import venv
import os
import os.path
import sys

join = os.path.join

if not "--virtualenv" in sys.argv:

    print("Creating new virtualenv with pip in {}".format(join(os.path.curdir, "env")))
    venv.create("env", with_pip=True)
    os.environ["PATH"] = os.path.realpath(join(os.path.curdir, "env", "bin")) + os.pathsep + os.environ["PATH"]

    print("Launching new python in virtualenv")
    os.execvp(os.path.basename(sys.executable), [os.path.basename(sys.executable), os.path.realpath(__file__), "--virtualenv"])
else:
    import pip
    print("Installing requirements")
    pip.main(["install", "-r", "requirements.txt"])

    from jupyter_core.command import main
    import os.path
    import re
    path_to_notebook = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Analysis.ipynb")
    sys.argv = [re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0]), 
            "nbextension", "enable", "--py", "widgetsnbextension", "--sys-prefix"]
    sys.exit(main())
