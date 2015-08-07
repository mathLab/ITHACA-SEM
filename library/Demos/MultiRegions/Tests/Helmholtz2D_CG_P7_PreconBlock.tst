<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 2D CG with P=7 and block preconditioner</description>
    <executable>Helmholtz2D</executable>
    <parameters>-v -I Preconditioner=Block Helmholtz2D_P7_Periodic.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_P7_Periodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">6.82372e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">9.43712e-07</value>
        </metric>
        <metric type="Precon" id="3">
            <value tolerance="2">14</value>
        </metric>
    </metrics>
</test>
