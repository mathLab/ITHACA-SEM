<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Hex Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Hex_Kovasnay_SubStep.xml</parameters>
    <files>
        <file description="Session File">Hex_Kovasnay_SubStep.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.048372</value>
            <value variable="v" tolerance="1e-6">0.00253099</value>
            <value variable="w" tolerance="1e-6">0.000796309</value>
            <value variable="p" tolerance="1e-6">0.00900028</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0616737</value>
            <value variable="v" tolerance="1e-6">0.00726816</value>
            <value variable="w" tolerance="1e-6">0.00562612</value>
            <value variable="p" tolerance="1e-6">0.0453188</value>
        </metric>
    </metrics>
</test>


