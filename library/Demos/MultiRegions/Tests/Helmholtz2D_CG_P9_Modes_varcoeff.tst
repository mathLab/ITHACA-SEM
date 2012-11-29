<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=9</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_CG_P9_Modes_varcoeff.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_CG_P9_Modes_varcoeff.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.00226302</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.00275842</value>
        </metric>
    </metrics>
</test>


