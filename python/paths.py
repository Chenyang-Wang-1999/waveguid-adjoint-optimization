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
import subprocess

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
    """拼出调用求解器的完整命令行（argv 列表，便于安全传递参数）。"""
    cmd = [MPIEXEC, "-n", "1", SOLVER]
    if extra_args:
        cmd += extra_args.split() if isinstance(extra_args, str) else list(extra_args)
    return cmd


def run_solver(extra_args="-pc_type lu"):
    """调用有限元求解器，返回退出码（0 表示成功）。

    返回的是真正的进程退出码，调用方应当检查它：求解器崩溃时如果继续去读
    上一次的结果文件，会拿到**旧的灵敏度**而不报错。
    """
    return subprocess.run(solver_command(extra_args)).returncode
