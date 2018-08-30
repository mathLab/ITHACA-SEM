<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_SubStep_2order.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_SubStep_2order.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.00511109</value>
            <value variable="v" tolerance="1e-6">0.00201302</value>
	    <value variable="p" tolerance="1e-6">0.0122811 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00959368</value>
            <value variable="v" tolerance="1e-6">0.00584268</value>
	    <value variable="p" tolerance="1e-6">0.0500805</value>
        </metric>
    </metrics>
</test>
