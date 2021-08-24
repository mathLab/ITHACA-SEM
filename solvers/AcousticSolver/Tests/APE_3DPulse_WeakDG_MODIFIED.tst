<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_3DPulse_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_3DPulse_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">6.48829</value>
            <value variable="u" tolerance="1e-12">0.00900547</value>
            <value variable="v" tolerance="1e-12">0.00900547</value>
            <value variable="w" tolerance="1e-12">0.00900547</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">4.28036</value>
            <value variable="u" tolerance="1e-5">0.0093</value>
            <value variable="v" tolerance="1e-5">0.0093</value>
            <value variable="w" tolerance="1e-5">0.0093</value>
        </metric>
    </metrics>
</test>
