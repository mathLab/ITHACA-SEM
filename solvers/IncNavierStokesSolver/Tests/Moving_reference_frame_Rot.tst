<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D flexible cylinder flow simulation using "MovingBody" module</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Moving_reference_frame_Rot.xml</parameters>
    <files>
        <file description="Session File">Moving_reference_frame_Rot.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0796359</value>
            <value variable="v" tolerance="1e-12">0.105918</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0.305458</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.243522</value>
            <value variable="v" tolerance="1e-12">0.499246</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">1.98065</value>
        </metric>
    </metrics>
</test>
