<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>APESolver</executable>
    <parameters>APE_1DPulseWall_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_1DPulseWall_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">000</value>
            <value variable="u" tolerance="1e-12">000</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">000</value>
            <value variable="u" tolerance="1e-12">000</value>
        </metric>
    </metrics>
</test>
