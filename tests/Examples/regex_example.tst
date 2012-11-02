<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Channel Flow P=2</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Channel_Flow_3modes.xml</parameters>
    <files>
        <file description="Session File">Channel_Flow_3modes.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>^L 2 error\s*(?:\(variable (\w+)\))?\s*:\s*([+-]?\d.+\d|0).*</regex>
            <matches>
                <match>
                    <field id="0">u</field>
                    <field id="1" tolerance="1e-12">5.07143e-16</field>
                </match>
                <match>
                    <field id="0">v</field>
                    <field id="1" tolerance="1e-12">0</field>
                </match>
                <match>
                    <field id="0">p</field>
                    <field id="1" tolerance="1e-8">7.35816e-15</field>
                </match>
            </matches>
        </metric>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.07143e-16</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-8">7.35816e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.66294e-15</value>
            <value variable="v" tolerance="1e-12">3.05946e-16</value>
            <value variable="p" tolerance="1e-8">5.75096e-14</value>
        </metric>
    </metrics>
</test>
