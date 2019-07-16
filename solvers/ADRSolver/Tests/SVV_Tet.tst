<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>SVV Tet P=6</description>
    <executable>ADRSolver</executable>
    <parameters>SVV_Tet.xml</parameters>
    <files>
        <file description="Session File">SVV_Tet.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="0">
            <value variable="u" tolerance="1e-05">0.000183629</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-05">0.00221343</value>
        </metric>
    </metrics>
</test>
