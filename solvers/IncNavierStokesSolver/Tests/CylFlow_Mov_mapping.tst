<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving cylinder, 2D, using Mappings for forced oscillations</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_Mov_mapping.xml</parameters>
    <files>
        <file description="Session File">CylFlow_Mov_mapping.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">21.4161</value>
            <value variable="v" tolerance="1e-9">0.911317</value>
	    <value variable="p" tolerance="1e-9">1.42587</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">1.87012</value>
            <value variable="v" tolerance="1e-9">0.92044</value>
	    <value variable="p" tolerance="1e-9">1.76508</value>
        </metric>
    </metrics>
</test>


