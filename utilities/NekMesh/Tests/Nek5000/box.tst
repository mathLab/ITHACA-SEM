<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nek5000 box</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list box.rea5000.gz box-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">box.rea5000.gz</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
