<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Channel Flow P=8, flowrate driven, parallel</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8_Flowrate.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File">ChanFlow_m8_Flowrate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.00299e-13</value>
            <value variable="v" tolerance="2e-12">1.84102e-13</value>
            <value variable="p" tolerance="1e-8">5.92475e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-12">5.93081e-13</value>
            <value variable="v" tolerance="4e-12">3.56319e-13</value>
            <value variable="p" tolerance="1e-8">9.17508e-10</value>
        </metric>
    </metrics>
</test>
