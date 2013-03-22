<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=7</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_CG_P7_Modes_curved_quad.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_CG_P7_Modes_curved_quad.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">6.31597e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">6.31597e-07</value>
        </metric>
    </metrics>
</test>


