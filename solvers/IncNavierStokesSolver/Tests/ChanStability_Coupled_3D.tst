<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Linear Stability (Arpack): Channel Flow</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanStability_Coupled_3D.xml</parameters>
    <files>
        <file description="Session File">ChanStability_Coupled_3D.xml</file>
	<file description="Session File">ChanStability_Coupled_3D.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">106244</value>
            <value variable="v" tolerance="1e-6">92.8801</value>
            <value variable="w" tolerance="1e-6">228.772</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0</value>
            <value variable="v" tolerance="1e-6">0</value>
            <value variable="w" tolerance="1e-6">0</value>
        </metric>
    </metrics>
</test>


