<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linear stability with coupled solver (Arpack LI): Channel Largest real Ev = (0.00223554,+/-0.249844i)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanStability_Coupled.xml</parameters>
    <files>
        <file description="Session File">ChanStability_Coupled.xml</file>
	<file description="Session File">ChanStability_Coupled.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="0.0001">2.58904</value>
            <value variable="v" tolerance="0.0001">0.00277094</value>
        </metric>
        <metric type="Linf" id="2">
		<value variable="u" tolerance="0.00001">1.00142</value>
            <value variable="v" tolerance="0.00001">0.00367426</value>
        </metric>
    </metrics>
</test>


