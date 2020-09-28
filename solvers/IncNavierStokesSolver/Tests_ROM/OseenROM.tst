<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=11</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>OseenROM.xml</parameters>
    <files>
        <file description="Session File">OseenROM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">16.9188</value>
            <value variable="v" tolerance="1e-12">1.07847</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6</value>
            <value variable="v" tolerance="1e-12">0.895333</value>
        </metric>
    </metrics>
</test>


