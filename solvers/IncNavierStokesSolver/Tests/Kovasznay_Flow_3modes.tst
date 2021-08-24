<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=2</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Kovasznay_Flow_3modes.xml</parameters>
    <files>
        <file description="Session File">Kovasznay_Flow_3modes.xml</file>
        <file description="Session File">Kovasznay_Flow_3modes.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-3">0.343606</value>
            <value variable="v" tolerance="1e-3">0.0255658</value>
            <value variable="p" tolerance="1e-3">0.347747</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-3">0.322039</value>
            <value variable="v" tolerance="1e-3">0.0243043</value>
            <value variable="p" tolerance="1e-3">0.519425</value>
        </metric>
    </metrics>
</test>


