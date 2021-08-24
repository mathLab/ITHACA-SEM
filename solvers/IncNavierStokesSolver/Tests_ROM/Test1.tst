<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=11</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test1.xml</parameters>
    <files>
        <file description="Session File">Test1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">16.919</value>
            <value variable="v" tolerance="1e-12">1.07889</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">6</value>
            <value variable="v" tolerance="1e-12">0.903853</value>
        </metric>
    </metrics>
</test>


