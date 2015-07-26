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
            <value variable="u" tolerance="0.01">14.5728</value>
            <value variable="v" tolerance="0.0001">0.0158214</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="0.01">6.91611</value>
            <value variable="v" tolerance="0.00001">0.0098673</value>
        </metric>
    </metrics>
</test>


