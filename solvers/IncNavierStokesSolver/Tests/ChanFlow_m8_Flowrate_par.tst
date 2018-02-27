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
            <value variable="u" tolerance="1e-12">7.66614e-11</value>
            <value variable="v" tolerance="1e-12">4.52553e-12</value>
            <value variable="p" tolerance="1e-8">7.0151e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.1156e-10</value>
            <value variable="v" tolerance="1e-12">1.17928e-11</value>
            <value variable="p" tolerance="1e-8">8.32526e-09</value>
        </metric>
    </metrics>
</test>
