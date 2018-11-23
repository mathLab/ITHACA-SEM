<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=400</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_2DPulseWall_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DPulseWall_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">6.7577</value>
            <value variable="u" tolerance="1e-7">0.014607</value>
            <value variable="v" tolerance="1e-7">0.00571716</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">13.659</value>
            <value variable="u" tolerance="1e-7">0.0280936</value>
            <value variable="v" tolerance="1e-7">0.0108992</value>
        </metric>
    </metrics>
</test>
