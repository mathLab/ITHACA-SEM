<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description> Phase sampling for time average filters</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-v PhaseSampling.xml</parameters>
    <files>
        <file description="Session File">PhaseSampling.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Sample: (\d+)\s*Phase: ([\d\.]+)\s*</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                    <field id="1">0.58</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
