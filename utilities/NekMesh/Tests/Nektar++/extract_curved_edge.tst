<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Extraction of curved surface</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list -m extract:surf=3 extract_curved_edge.xml extract_curved_edge-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">extract_curved_edge.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>.*Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
