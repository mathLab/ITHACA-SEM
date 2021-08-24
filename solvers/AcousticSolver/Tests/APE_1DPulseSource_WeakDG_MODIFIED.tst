<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_1DPulseSource_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_1DPulseSource_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">22.2828</value>
            <value variable="u" tolerance="1e-12">0.0748752</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">27.7849</value>
            <value variable="u" tolerance="1e-12">0.104506</value>
        </metric>
    </metrics>
</test>
