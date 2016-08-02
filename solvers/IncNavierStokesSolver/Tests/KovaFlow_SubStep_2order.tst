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
            <value variable="u" tolerance="1e-6">0.00504022</value>
            <value variable="v" tolerance="1e-6">0.00200113</value>
	    <value variable="p" tolerance="1e-6">0.0123557 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00959985</value>
            <value variable="v" tolerance="1e-6">0.00582256</value>
	    <value variable="p" tolerance="1e-6">0.0502862</value>
        </metric>
    </metrics>
</test>
