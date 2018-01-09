<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=10</description>
    <executable>APESolver</executable>
    <parameters>APE_2DVariableC_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DVariableC_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">52.0346</value>
            <value variable="u" tolerance="1e-7">0.228765</value>
            <value variable="v" tolerance="1e-7">0.0598651</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">60.1152</value>
            <value variable="u" tolerance="1e-7">0.246703</value>
            <value variable="v" tolerance="1e-7">0.109342</value>
        </metric>
    </metrics>
</test>
