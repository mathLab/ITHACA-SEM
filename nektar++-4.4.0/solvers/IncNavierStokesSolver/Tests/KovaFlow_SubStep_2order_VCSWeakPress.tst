<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_SubStep_2order_VCSWeakPress.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_SubStep_2order_VCSWeakPress.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.00494118</value>
            <value variable="v" tolerance="1e-6">0.00193291</value>
	    <value variable="p" tolerance="1e-6">0.0130327 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0093622</value>
            <value variable="v" tolerance="1e-6">0.00584567</value>
	    <value variable="p" tolerance="1e-6">0.0513235</value>
        </metric>
    </metrics>
</test>
