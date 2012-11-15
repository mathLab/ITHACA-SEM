<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=2</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_Flow_3modes.xml</parameters>
    <files>
        <file description="Session File">Channel_Flow_3modes.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.20394e-16</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">8.50369e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.7462e-15</value>
            <value variable="v" tolerance="1e-12">3.04438e-16</value>
            <value variable="p" tolerance="1e-12">6.23945e-14</value>
        </metric>
    </metrics>
</test>


