<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=7 XXT Full</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_CG_P7_Modes_AllBCs_xxt_full.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0.00888081</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.0101781</value>
        </metric>
    </metrics>
</test>


