"""生成/重建一个算例：gmsh 网格 + design.h5。

用法::

    python3 ../python/make_mesh.py

默认沿用 ``data/DOF_list.csv`` 最后一轮的设计变量（断点续跑时最常用）；
没有历史记录时用初始的直线设计。
"""

import os

import numpy as np

import design
import optimizer
from paths import chdir_data


def main(name="cross_waveguide", n_ctrl=10):
    chdir_data()

    dofs = optimizer.last_dofs()
    if dofs is not None:
        source = "%s 的最后一轮" % optimizer.DOF_CSV
    else:
        geo, layout = design.Geometry(), design.Layout(n_ctrl)
        dofs = design.initial_dofs(geo, layout)
        source = "初始直线设计"

    print("设计变量来源：%s（%d 个自由度）" % (source, len(dofs)))
    case = design.make_case(dofs, name=name)
    print("网格: %s" % case.mesh_file)
    print("设计: %s" % case.design_file)
    if case.ids:
        print("physical id: %s" % case.ids)
    return case


if __name__ == "__main__":
    main()
