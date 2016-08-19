<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Channel Flow P=8, flowrate driven</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m8_Flowrate.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m8_Flowrate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.81544e-16</value>
            <value variable="v" tolerance="1e-12">2.00905e-16</value>
            <value variable="p" tolerance="1e-8">9.32789e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.03806e-14</value>
            <value variable="v" tolerance="1e-12">1.02362e-15</value>
            <value variable="p" tolerance="1e-8">3.4639e-14</value>
        </metric>
    </metrics>
</test>
