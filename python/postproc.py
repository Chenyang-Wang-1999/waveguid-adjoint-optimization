"""后处理：读场、可视化。

分两个层次，别混用：

* :func:`read_field` / :func:`plot` —— **只读已有结果**，只需要文件名，
  不需要 :class:`design.Case`。日常画图用这个。
* :func:`post_process` / :func:`visualization` —— 需要 ``save_field=1`` 的场，
  而优化循环里是 ``save_field=0``，所以它们会**重跑一次求解器**，因此要 Case。

例::

    import postproc
    postproc.plot("out_cross_waveguide.msh", "result.h5")

    # 或者：从 design.h5 还原算例，重跑一次带场的正问题再画
    import design
    postproc.visualization(design.Case.load())
"""

import numpy as np
import scipy.io as spio

import dataio
import solver
import paths
from design import Case

DEFAULT_RESULT = "result.h5"


def read_field(result_file=DEFAULT_RESULT):
    """从已经算好的结果里读场，不重跑求解器。"""
    res = dataio.read_result(result_file, want_field=True)
    if "Ez" not in res:
        raise RuntimeError(
            "%s 里没有 /fields/Ez；求解时需要把 infile.pot 的 save_field 设为 1"
            % result_file)
    return res["Ez"]


def plot(mesh_file, result_file=DEFAULT_RESULT, scalars="Ez", Ez=None):
    """把场画出来（只读文件，不跑求解器）。"""
    import pyvista

    if Ez is None:
        Ez = read_field(result_file)
    mesh = pyvista.read(mesh_file)
    mesh["E_norm"] = np.abs(Ez)
    mesh["Ez"] = Ez.real
    return mesh.plot(scalars=scalars)


def post_process(case, out_mat=None):
    """重跑一次正问题（``save_field=1``）并把场导出成 .mat。

    只有这个（和下面的 :func:`visualization`）需要 ``case``，因为要真的
    再解一遍才能拿到场。
    """
    res = solver.run(case, save_field=True)
    if "Ez" not in res:
        raise RuntimeError("结果里没有场；确认求解器成功写出了 /fields")
    spio.savemat(out_mat or case.name + "_post.mat", {"Ez": res["Ez"]})
    return res["Ez"]


def visualization(case, scalars="Ez", out_mat=None):
    """重跑一次并作图（旧 ``visualization()`` 的行为）。"""
    Ez = post_process(case, out_mat=out_mat)
    return plot(case.out_mesh_file, case.result_file, scalars=scalars, Ez=Ez)


if __name__ == "__main__":
    paths.chdir_data()
    cs = Case.load("design.h5")
    print(cs)
    visualization(cs)