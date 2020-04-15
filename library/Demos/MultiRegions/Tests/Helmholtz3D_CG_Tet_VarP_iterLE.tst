<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 3D CG for Tet with variable P</description>
    <executable>Helmholtz3D</executable>
    <parameters>-I GlobalSysSoln=IterativeStaticCond -I Preconditioner=LowEnergyBlock Helmholtz3D_Tet_VarP.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Tet_VarP.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.00416454</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.0313733</value>
        </metric>
    </metrics>
</test>


