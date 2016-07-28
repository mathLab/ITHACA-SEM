<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Channel Flow P=8, flowrate driven, parallel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">ChanFlow_m8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">6.13359e-09</value>
            <value variable="v" tolerance="1e-12">5.74846e-09</value>
            <value variable="p" tolerance="1e-8">1.37568e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.02619e-07</value>
            <value variable="v" tolerance="1e-12">2.55204e-08</value>
            <value variable="p" tolerance="1e-8">5.02791e-07</value>
        </metric>
    </metrics>
</test>
