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
            <value variable="u" tolerance="1e-12">9.73822e-13</value>
            <value variable="v" tolerance="1e-12">1.29647e-12</value>
            <value variable="p" tolerance="1e-8">9.23756e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.33813e-12</value>
            <value variable="v" tolerance="1e-12">2.31607e-12</value>
            <value variable="p" tolerance="1e-8">1.23396e-09</value>
        </metric>
    </metrics>
</test>
