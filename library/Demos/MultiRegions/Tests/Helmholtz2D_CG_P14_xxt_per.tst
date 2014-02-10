<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=14 xxt periodic</description>
    <executable>Helmholtz2D</executable>
    <parameters>Helmholtz2D_CG_P14_xxt_per.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_CG_P14_xxt_per.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">6.22012e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">8.34899e-12</value>
        </metric>
    </metrics>
</test>


