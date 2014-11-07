<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>APESolver</executable>
    <parameters>APE_3DPulseWall_FRDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_3DPulseWall_FRDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">6.32197</value>
            <value variable="u" tolerance="1e-7">0.00919153</value>
            <value variable="v" tolerance="1e-7">0.00865487</value>
            <value variable="w" tolerance="1e-7">0.00865487</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">11.9004</value>
            <value variable="u" tolerance="1e-7">0.0213414</value>
            <value variable="v" tolerance="1e-7">0.0212787</value>
            <value variable="w" tolerance="1e-7">0.0212788</value>
        </metric>
    </metrics>
</test>
