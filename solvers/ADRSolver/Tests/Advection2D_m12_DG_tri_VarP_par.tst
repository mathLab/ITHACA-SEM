<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady DG advection, tri, order 4, P=Variable</description>
    <executable>ADRSolver</executable>
    <parameters>--use-scotch Advection2D_m12_DG_tri_VarP.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">Advection2D_m12_DG_tri_VarP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.03599e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.24106e-05</value>
        </metric>
    </metrics>
</test>


