"""Functions for trees."""

from string import whitespace

escaped_tokens = set('\\')
structure_tokens = set('(),:;')
whitespace_tokens = set(whitespace)


class TreeNode:
    """
    A class representing nodes on a tree.

    Parameters
    ----------
    name: str
        Name of node.
    children: list of TreeNodes
        Children of node.
    parent: TreeNode
        Parent of node.
    length: Length of node.
    """

    def __init__(self, name=None, children=None, parent=None, length=None):
        if children is None:
            children = []
        self.name = name
        self.children = children
        self.parent = parent
        self.length = length
        for child in self.children:
            child.parent = self

    def __repr__(self):
        class_name = type(self).__name__
        node_name = self.name if self.name is not None else 'unnamed'
        children_count = len(self.children)
        repr_string = f'<{class_name}, name: {node_name}, children count: {children_count}>'
        return repr_string

    def __str__(self):
        return self.to_newick()

    def traverse(self, include_self=True, order='pre'):
        """
        Returns an iterator over nodes.

        Parameters
        ----------
        include_self : bool
            If True, include calling node in iterator.
        order : str
            Type of traversal. 'pre', 'post', or 'level' options are supported.

        Returns
        -------
        iterator over TreeNodes
        """
        if order == 'pre':
            return self._preorder(include_self=include_self)
        elif order == 'post':
            return self._postorder(include_self=include_self)
        elif order == 'level':
            return self._levelorder(include_self=include_self)
        else:
            raise ValueError('order is not one of pre, post, or level')

    def _levelorder(self, include_self=True):
        stack = []
        if include_self:
            stack.append(self)
        else:
            stack.extend(self.children)
        while stack:
            current_node = stack.pop(0)
            yield current_node
            stack.extend(current_node.children)

    def _postorder(self, include_self=True):
        stack = []
        if include_self:
            stack.append(self)
        else:
            stack.extend(self.children)
        while stack:
            current_node = stack.pop()
            if stack and current_node is stack[-1]:
                yield stack.pop()
            else:
                stack.extend([current_node, current_node])
                stack.extend(current_node.children[::-1])

    def _preorder(self, include_self=True):
        stack = []
        if include_self:
            stack.append(self)
        else:
            stack.extend(self.children[::-1])
        while stack:
            current_node = stack.pop()
            yield current_node
            stack.extend(current_node.children[::-1])

    def is_root(self):
        return self.parent is None

    def is_tip(self):
        return not bool(self.children)

    def tips(self, include_self=True, order='pre'):
        for node in self.traverse(include_self=include_self, order=order):
            if node.is_tip():
                yield node

    def non_tips(self, include_self=True, order='pre'):
        for node in self.traverse(include_self=include_self, order=order):
            if not node.is_tip():
                yield node

    def copy(self):
        copies = []
        for node in self.traverse(order='post'):
            if node.is_tip():
                cls = type(node)
                copy = cls(name=node.name, length=node.length)
                copies.append(copy)
            else:
                cls = type(node)
                children_copies = [copies.pop() for _ in node.children][
                    ::-1
                ]  # Reverse b/c post order reverses children
                copy = cls(name=node.name, children=children_copies, length=node.length)
                copies.append(copy)
        copy = copies.pop()
        return copy

    @staticmethod
    def from_newick(nw_string):
        nw_tokens = _tokenize_newick(nw_string)
        tree = _parse_newick_tokens(nw_tokens)
        return tree

    def to_newick(self):
        nw_strings = []
        for node in self.traverse(order='post'):
            name = node.name if node.name is not None else ''
            length = f':{node.length}' if node.length else ''
            if node.is_tip():
                nw_strings.append(f'{name}{length}')
            else:
                children_newicks = [nw_strings.pop() for _ in node.children]
                nw_string = '(' + ', '.join(children_newicks) + f'){name}{length}'
                nw_strings.append(nw_string)
        nw_string = nw_strings.pop()
        return nw_string + ';'


def read_newick(path):
    """
    Read Newick file at path and return tree.

    Parameters
    ----------
    path : str
        Path to Newick file

    Returns
    -------
    tree : TreeNode
    """
    with open(path) as file:
        return parse_newick(file)


def parse_newick(file):
    """
    Read Newick from stream and return tree.

    Parameters
    ----------
    file : file object
        File object of Newick text stream.

    Returns
    -------
    tree : TreeNode
    """
    nw_string = file.read().rstrip()
    tree = TreeNode.from_newick(nw_string)
    return tree


def write_newick(tree, path):
    """
    Write tree to file in Newick format.

    Parameters
    ----------
    tree : TreeNode
    path : str
    """
    nw_string = tree.to_newick()
    with open(path, 'w') as file:
        file.write(nw_string)


def _tokenize_newick(nw_string):
    is_escaped = False
    tokens = []
    buffer = []  # Buffer for label or length

    i = 0
    while i < len(nw_string):
        sym = nw_string[i]
        if not is_escaped and sym in escaped_tokens:
            is_escaped = True
            i += 1
        elif not is_escaped and sym in structure_tokens:
            if buffer:
                tokens.append(''.join(buffer))
                buffer = []
            tokens.append(sym)
            i += 1
        elif not is_escaped and sym in whitespace_tokens:
            if buffer:
                tokens.append(''.join(buffer))
                buffer = []
            i += 1
        else:  # Otherwise sym is part of label or length
            buffer.append(sym)
            is_escaped = False
            i += 1

    return tokens


def _parse_newick_tokens(nw_tokens):
    current_node = TreeNode()
    current_parent = None

    i = 0
    while i < len(nw_tokens):
        token = nw_tokens[i]
        if token == '(':
            current_parent = current_node
            current_node = TreeNode(parent=current_parent)
            current_parent.children.append(current_node)
            i += 1
        elif token == ')':
            current_node = current_node.parent
            current_parent = current_node.parent
            i += 1
        elif token == ',':
            current_node = TreeNode(parent=current_parent)
            current_parent.children.append(current_node)
            i += 1
        elif token == ':':
            length_token = nw_tokens[i + 1]
            current_node.length = float(length_token)
            i += 2
        elif token == ';':
            i += 1
        else:
            current_node.name = token
            i += 1

    return current_node
