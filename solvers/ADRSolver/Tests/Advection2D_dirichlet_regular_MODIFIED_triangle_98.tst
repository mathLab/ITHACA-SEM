<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady WeakDG advection MODIFIED, P=3 Dirichlet bcs, regular triangular elements using AVX BwdTrans</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_dirichlet_regular_MODIFIED_triangle_98.xml</parameters>
    <files>
        <file description="Session File">Advection2D_dirichlet_regular_MODIFIED_triangle_98.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7"> 0.0014953 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8"> 0.00385391 </value>
        </metric>
    </metrics>
</test>
