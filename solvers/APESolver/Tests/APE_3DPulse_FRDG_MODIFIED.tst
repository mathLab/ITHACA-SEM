<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>APESolver</executable>
    <parameters>APE_3DPulse_FRDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_3DPulse_FRDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">6.48826</value>
            <value variable="u" tolerance="1e-7">0.00900561</value>
            <value variable="v" tolerance="1e-7">0.00900546</value>
            <value variable="w" tolerance="1e-7">0.00900547</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">4.28815</value>
            <value variable="u" tolerance="1e-7">0.00930734</value>
            <value variable="v" tolerance="1e-7">0.00930633</value>
            <value variable="w" tolerance="1e-7">0.00930837</value>
        </metric>
    </metrics>
</test>
