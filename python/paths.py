"""项目路径与外部命令的统一入口。

约定：所有脚本都在 ``data/`` 目录下工作。有限元求解器从*当前目录*读 ``infile.pot``，
并把 ``out_*.msh`` / ``*_result.mat`` 写在同一个目录，所以脚本开头调用一次
:func:`chdir_data` 即可，其余相对路径写法保持不变。

路径都可以用环境变量覆盖，便于在别的机器上跑：

===========================  ==============================
环境变量                     含义
===========================  ==============================
``WAVEGUIDE_SOLVER``         有限元求解器可执行文件的完整路径
``MPIEXEC``                  mpiexec 的路径（默认从 PATH 里找）
``WAVEGUIDE_DATA_DIR``       data 目录的位置
===========================  ==============================
"""

import os
import shutil

PYTHON_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(PYTHON_DIR)

DATA_DIR = os.environ.get("WAVEGUIDE_DATA_DIR", os.path.join(PROJECT_ROOT, "data"))

SOLVER = os.environ.get(
    "WAVEGUIDE_SOLVER",
    os.path.join(PROJECT_ROOT, "waveguide_optimization", "waveguide_optimization"),
)

MPIEXEC = os.environ.get("MPIEXEC") or shutil.which("mpiexec") or "mpiexec"


def chdir_data():
    """切到 data/ 目录（没有就建一个）。"""
    os.makedirs(DATA_DIR, exist_ok=True)
    os.chdir(DATA_DIR)


def solver_command(extra_args="-pc_type lu"):
    """拼出调用求解器的完整命令行。"""
    cmd = "%s -n 1 %s" % (MPIEXEC, SOLVER)
    if extra_args:
        cmd += " " + extra_args
    return cmd


def run_solver(extra_args="-pc_type lu"):
    """调用有限元求解器，返回 shell 退出码。"""
    return os.system(solver_command(extra_args))
