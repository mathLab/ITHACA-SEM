<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Extrusion of 2D mesh</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list -m extrude:layers=1:length=0.1 extrude.xml
        extrude-out.xml:xml:test:stats</parameters>
    <files>
        <file description="Input File">extrude.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians\s*:\s*(\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="2">
            <regex>^Node count\s*:\s*(\d+)</regex>
            <matches>
                <match>
                    <field id="0">12</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="3">
            <regex>^Number of composites\s*:\s*(\d+)</regex>
            <matches>
                <match>
                    <field id="0">8</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="4">
            <regex>^Lower mesh extent\s*:\s*(-?\d+\.?\d?+)\s*(-?\d+\.?\d?+)\s*
                (-?\d+\.?\d?+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                    <field id="1">-0.1</field>
                    <field id="2">0</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="5">
            <regex>^Upper mesh extent\s*:\s*(-?\d+\.?\d?+)\s*(-?\d+\.?\d?+)\s*
                (-?\d+\.?\d?+)</regex>
            <matches>
                <match>
                    <field id="0">1</field>
                    <field id="1">1</field>
                    <field id="2">0.1</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>