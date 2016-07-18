<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=400</description>
    <executable>APESolver</executable>
    <parameters>APE_2DPulseWall_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DPulseWall_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">7.56584</value>
            <value variable="u" tolerance="1e-7">0.0174419</value>
            <value variable="v" tolerance="1e-7">0.00907637</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">18.7643</value>
            <value variable="u" tolerance="1e-7">0.0322549</value>
            <value variable="v" tolerance="1e-7">0.0460461</value>
        </metric>
    </metrics>
</test>
