<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Helmholtz 2D CG with P=5 with disconnected domain</description>
    <executable>Helmholtz2D</executable>
    <parameters>-v Helmholtz2D_Disconnected.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_Disconnected.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">3.24525e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.35327e-07</value>
        </metric>
    </metrics>
</test>
