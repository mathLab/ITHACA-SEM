<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 1D HDG with P=8</description>
    <executable>HDGHelmholtz1D</executable>
    <parameters>helmholtz1D_8modes.xml</parameters>
    <files>
        <file description="Session File">helmholtz1D_8modes.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">5.81155e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">0.000213838</value>
        </metric>
    </metrics>
</test>


