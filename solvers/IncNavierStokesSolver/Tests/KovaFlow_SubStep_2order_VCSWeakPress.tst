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
            <value variable="u" tolerance="1e-6">0.00494845</value>
            <value variable="v" tolerance="1e-6">0.00193291</value>
	    <value variable="p" tolerance="1e-6">0.0130356 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00936783</value>
            <value variable="v" tolerance="1e-6">0.005867</value>
	    <value variable="p" tolerance="1e-6">0.051281</value>
        </metric>
    </metrics>
</test>
