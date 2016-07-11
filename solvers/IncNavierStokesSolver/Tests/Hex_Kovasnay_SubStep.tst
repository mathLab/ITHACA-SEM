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
            <value variable="u" tolerance="1e-6">0.0525209 </value>
            <value variable="v" tolerance="1e-6">0.00347838</value>
            <value variable="w" tolerance="1e-6">0.00339679</value>
            <value variable="p" tolerance="1e-6">0.0136489 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0637047 </value>
            <value variable="v" tolerance="1e-6">0.00776541</value>
            <value variable="w" tolerance="1e-6">0.00629253</value>
            <value variable="p" tolerance="1e-6">0.0470333 </value>
        </metric>
    </metrics>
</test>
