<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (Arpack): Channel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanStability_Coupled.xml</parameters>
    <files>
        <file description="Session File">ChanStability_Coupled.xml</file>
	<file description="Session File">ChanStability_Coupled.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="100">2.58948</value>
            <value variable="v" tolerance="100">0.000352934</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="100">1.01298</value>
            <value variable="v" tolerance="100">0.000295924</value>
        </metric>
    </metrics>
</test>


