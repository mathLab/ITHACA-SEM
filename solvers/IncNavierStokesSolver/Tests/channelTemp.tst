<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=4 with temperature field</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>channelTemp.xml</parameters>
    <files>
        <file description="Session File">channelTemp.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">8.82097e-07</value>
            <value variable="v" tolerance="1e-10">1.56407e-07</value>
	    <value variable="c1"tolerance="1e-4"> 1.02325</value>
            <value variable="p" tolerance="1e-3"> 0.048 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.12525e-05</value>
            <value variable="v" tolerance="1e-8">1.31727e-06</value>
	    <value variable="c1"tolerance="1e-4">0.947434</value>
	    <value variable="p" tolerance="1e-4">0.0240083</value>
        </metric>
    </metrics>
</test>

