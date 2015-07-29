<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 2D CG with P=7 and block preconditioner</description>
    <executable>Helmholtz2D</executable>
    <parameters>-v Helmholtz2D_CG_PreconBlock_P7.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_CG_PreconBlock_P7.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">6.82374e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">9.43919e-07</value>
        </metric>
        <metric type="Precon" id="3">
            <value tolerance="2">17</value>
        </metric>
    </metrics>
</test>
