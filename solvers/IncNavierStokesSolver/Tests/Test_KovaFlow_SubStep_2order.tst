<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasnay solution using sub-stepping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_KovaFlow_Substep_2order.xml</parameters>
    <files>
        <file description="Session File">Test_KovaFlow_SubStep_2order.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.00528173</value>
            <value variable="v" tolerance="1e-12">0.0020437</value>
	    <value variable="p" tolerance="1e-12">0.0123206</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.00999388</value>
            <value variable="v" tolerance="1e-12">0.00586634</value>
	    <value variable="p" tolerance="1e-12">0.0501208</value>
        </metric>
    </metrics>
</test>