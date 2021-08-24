<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Rayleigh Bernard flow at P=4 </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>rbc.xml</parameters>
    <files>
        <file description="Session File">rbc.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">2.31672e-10</value>
            <value variable="v" tolerance="1e-11">7.86654e-11</value>
	    <value variable="T"tolerance="1e-4"> 0.820569</value>
            <value variable="p" tolerance="1e-3"> 1128.21 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">3.78998e-10</value>
            <value variable="v" tolerance="1e-10">1.32329e-10</value>
	    <value variable="T"tolerance="1e-2">1.0</value>
	    <value variable="p" tolerance="1e-2">1775</value>
        </metric>
    </metrics>
</test>

