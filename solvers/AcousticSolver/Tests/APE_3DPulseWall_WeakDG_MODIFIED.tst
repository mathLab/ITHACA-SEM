<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=300</description>
    <executable>AcousticSolver</executable>
    <parameters>APE_3DPulseWall_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_3DPulseWall_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">6.32197</value>
            <value variable="u" tolerance="1e-12">0.00919153</value>
            <value variable="v" tolerance="1e-12">0.00865486</value>
            <value variable="w" tolerance="1e-12">0.00865486</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">11.9003</value>
            <value variable="u" tolerance="1e-12">0.0213413</value>
            <value variable="v" tolerance="1e-12">0.0212786</value>
            <value variable="w" tolerance="1e-12">0.0212786</value>
        </metric>
    </metrics>
</test>
